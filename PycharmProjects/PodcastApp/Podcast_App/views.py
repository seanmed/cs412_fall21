from django.shortcuts import render
from django.http import HttpResponse
from .forms import SearchSpotifyAPIForm
# Create your views here.

def contact(request):
    if request.method == 'POST':
        form = SearchSpotifyAPIForm(request.POST)
        if form.is_valid():
            twitter_username= form.cleaned_data['twitter_username']
            print(twitter_username)


    form= SearchSpotifyAPIForm()
    return render(request, 'form.html', {'form': form})
